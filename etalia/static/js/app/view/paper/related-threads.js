define([
    'app/app',
    'text!app/templates/paper/related-threads.hbs',
    'app/view/detail',
    'app/view/thread/thumb',
    'app/view/thread/detail'
], function (App, template, Detail) {

    App.View.Paper = App.View.Paper || {};

    App.View.Paper.RelatedThreads = App.Backbone.View.extend({
        tagName: 'div',
        attributes: {
            'class': 'neighbors-thumbs'
        },

        template: App.Handlebars.compile(template),

        paperId: null,
        buttons: null,
        returnCallback: null,

        thumbPrefix: 'paper-related-threads-thumb-',
        activeTimespan: null,
        listView: null,

        events: {
            "click .neighbors-timespan-selector a": "onTimespanSelectorClick"
        },

        initialize: function (options) {
            this.paperId = parseInt(options.paper_id);
            if (!this.paperId) {
                throw 'options.paper_id is mandatory';
            }
            if (!options.buttons) {
                throw 'options.buttons is mandatory';
            }
            this.buttons = options.buttons;
            this.returnCallback = options.return_callback;
            if (!(typeof this.returnCallback == 'function')) {
                throw 'options.return_callback is mandatory';
            }

            this.collection = new App.Model.Threads();
            this.collection.url = App.config.api_root + '/library/my-papers/' + this.paperId + '/related-threads';

            this.collection.on("add", this.onCollectionAdd, this);
            this.collection.on("remove", this.onCollectionRemove, this);

            this.listenTo(this, "model:detail", this.openDetail);

            // Collection
            //this.listenTo(this.collection, "add", this.onCollectionAdd);
            //this.listenTo(this.collection, "remove", this.onCollectionRemove);
        },

        openDetail: function(model) {
            var options = {
                model: model,
                buttons: this.buttons
            };
            var detailModel = new App.Model.Detail({
                view: new App.View.Thread.Detail(options)
            });
            detailModel.setCenterButton({
                icon: 'close',
                title: 'Back to previous paper',
                callback: this.returnCallback
            });

            model
                .fetch({data: {view: 'nested'}})
                .done(function() {
                    Detail.setModel(detailModel);
                });
        },

        render: function () {
            App.log('PaperRelatedThreadsView::render');

            this.$el.html(this.template({}));

            // Thumbs list
            this.listView = new App.View.List.create({}, {
                $target: this.$('div[data-related-threads-list-placeholder]')
            });
            this.pushSubView(this.listView);

            // Initial timespan
            this.activeTimespan = parseInt(this.$('.neighbors-timespan-selector li.active a').data('value'));

            this.load();

            return this;
        },

        load: function() {
            this.collection.reset();
            this.clearList();

            this.collection.setQuery({
                'time_span': this.activeTimespan
            });

            this.collection.fetch();
        },

        onTimespanSelectorClick: function(e) {
            e.preventDefault();

            this.$('.neighbors-timespan-selector li').removeClass('active');

            var $link = $(e.target).closest('a');
            this.activeTimespan = parseInt($link.data('value'));

            $link.closest('li').addClass('active');

            this.load();
        },

        onCollectionAdd: function (model) {
            //App.log('PaperRelatedThreadsView::onCollectionAdd');

            // Render the thumb
            var thumbView = new App.View.Thread.Thumb({
                id: this.thumbPrefix + model.get('id'),
                model: model,
                list: this
            });

            this.listView.addThumbView(thumbView);
        },

        onCollectionRemove: function (model) {
            //App.log('PaperRelatedThreadsView::onCollectionRemove');

            this.listView.removeThumbById(this.thumbPrefix + model.get('id'));
        },

        clearList: function() {
            this.listView.clear();
        }
    });


    App.View.Paper.RelatedThreads.create = function(options, createOptions) {
        options = options || {};

        var view = new App.View.Paper.RelatedThreads(options);
        if (createOptions) {
            App.View.create(view, createOptions);
        }

        return view;
    };

    return App.View.Paper.RelatedThreads;
});
