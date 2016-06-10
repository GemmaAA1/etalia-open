define([
    'app',
    'text!app/templates/thread/neighbors.hbs'
], function (App, template) {

    App.View.Thread = App.View.Thread || {};

    App.View.Thread.Neighbors = App.Backbone.View.extend({
        tagName: 'div',
        attributes: {
            'class': 'neighbors-thumbs'
        },

        template: App.Handlebars.compile(template),

        threadId: null,
        activeTimespan: null,

        thumbPrefix: 'thread-neighbors-thumb-',
        listView: null,

        events: {
            "click .neighbors-timespan-selector a": "onTimespanSelectorClick"
        },

        initialize: function (options) {
            this.threadId = parseInt(options.thread_id);
            if (!this.threadId) {
                throw 'options.thread_id is mandatory';
            }

            this.collection = new App.Model.Threads();
            this.collection.url = App.config.api_root + '/thread/threads/' + this.threadId + '/neighbors';

            //this.collection.on("add", this.onCollectionAdd, this);
            //this.collection.on("remove", this.onCollectionRemove, this);

            // Collection
            this.listenTo(this.collection, "add", this.onCollectionAdd);
            this.listenTo(this.collection, "remove", this.onCollectionRemove);
        },

        onTimespanSelectorClick: function(e) {
            e.preventDefault();

            this.$('.neighbors-timespan-selector li').removeClass('active');

            var $link = $(e.target).closest('a');
            this.activeTimespan = parseInt($link.data('value'));

            $link.closest('li').addClass('active');

            this.load();
        },

        render: function () {
            App.log('ThreadNeighborsView::render');

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

        onCollectionAdd: function (model) {
            App.log('ThreadNeighborsView::onCollectionAdd');

            // Render the thumb
            var thumbView = new App.View.Thread.Thumb({
                id: this.thumbPrefix + model.get('id'),
                model: model,
                list: this
            });

            this.listView.addThumbView(thumbView);
        },

        onCollectionRemove: function (model) {
            App.log('ThreadNeighborsView::onCollectionRemove');

            this.listView.removeThumbById(this.thumbPrefix + model.get('id'));
        },

        clearList: function() {
            this.listView.clear();
        }
    });


    App.View.Thread.Neighbors.create = function(options, createOptions) {
        options = options || {};

        var view = new App.View.Thread.Neighbors(options);
        if (createOptions) {
            App.View.create(view, createOptions);
        }

        return view;
    };

    return App.View.Thread.Neighbors;
});