define([
    'app',
    'text!app/templates/paper/neighbors.hbs',
    'app/view/detail'
], function (App, template, Detail) {

    App.View.Paper = App.View.Paper || {};

    App.View.Paper.Neighbors = App.Backbone.View.extend({
        tagName: 'div',
        attributes: {
            'class': 'neighbors-thumbs'
        },

        template: App.Handlebars.compile(template),

        paperId: null,
        buttons: null,
        returnCallback: null,

        thumbPrefix: 'paper-neighbors-thumb-',
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

            this.collection = new App.Model.Papers();
            this.collection.url = App.config.api_root + '/library/papers/' + this.paperId + '/neighbors';

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
                view: new App.View.Paper.Detail(options)
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
            App.log('PaperNeighborsView::render');

            this.$el.html(this.template({}));

            // Thumbs list
            this.listView = new App.View.List.create({}, {
                $target: this.$('div[data-neighbors-list-placeholder]')
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
                'time-span': this.activeTimespan
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
            //App.log('PaperNeighborsView::onCollectionAdd');

            // Render the thumb
            var thumbView = new App.View.Paper.Thumb({
                id: this.thumbPrefix + model.get('id'),
                model: model,
                list: this
            });

            this.listView.addThumbView(thumbView);
        },

        onCollectionRemove: function (model) {
            //App.log('PaperNeighborsView::onCollectionRemove');

            this.listView.removeThumbById(this.thumbPrefix + model.get('id'));
        },

        clearList: function() {
            this.listView.clear();
        }
    });


    App.View.Paper.Neighbors.create = function(options, createOptions) {
        options = options || {};

        var view = new App.View.Paper.Neighbors(options);
        if (createOptions) {
            App.View.create(view, createOptions);
        }

        return view;
    };

    return App.View.Paper.Neighbors;
});