define([
    'app',
    'text!app/templates/thread/list.html',
    'app/view/detail',
    'app/view/thread/detail',
    'app/view/thread/thumb',
    'app/collection/thread/thread'
], function (App, template) {

    var defaults = {
        query: {
            view: 'nested'
        }
    };

    App.View.Thread = App.View.Thread || {};

    return App.View.Thread.List = App.Backbone.View.extend({

        //tagName: 'div',
        //className: 'thumb-list',
        template: App.Handlebars.compile(template),

        detailView: null,
        $thumbList: null,

        events: {
            "click #thread-next-page": "onNextPageClick"
        },

        initialize: function (options) {
            App.defaults(options, defaults);

            this.collection = new App.Collection.Threads({
                query: options.query
            });

            this.listenTo(this.collection.fullCollection, "reset", this.onCollectionReset);
            this.listenTo(this.collection.fullCollection, "add", this.onCollectionAdd);
            this.listenTo(this.collection.fullCollection, "remove", this.onCollectionRemove);

            this.listenTo(this, "model:remove", this.onModelRemove);
            this.listenTo(this, "model:detail", this.openDetail);
        },

        openDetail: function (model) {
            App.log('ThreadListView::onModelDetail');

            var detailModel = new App.Model.Detail({
                list: this,
                view: new App.View.Thread.Detail({
                    model: model
                })
            });

            if (this.detailView) {
                this.detailView.model = detailModel;
            } else {
                this.detailView = new App.View.Detail({
                    model: detailModel
                })
            }

            this.detailView.render();
        },

        onModelRemove: function (model) {
            App.log('ThreadListView::onModelRemove');

            this.collection.remove(model);
            this.collection.fullCollection.remove(model);
        },

        onCollectionAdd: function (model) {
            App.log('ThreadListView::onCollectionAdd');

            // Render the thumb
            var thumbView = new App.View.Thread.Thumb({
                    id: 'thread-thumb-' + model.get('id'),
                    model: model,
                    list: this
                });

            // Append the thumb
            thumbView
                .render()
                .$el.appendTo(this.$thumbList);
        },

        onCollectionRemove: function (model) {
            App.log('ThreadListView::onCollectionRemove');
            this.$('#thread-thumb-' + model.get('id')).remove();
        },

        onCollectionReset: function () {
            if (this.collection.fullCollection.length == 0) {
                this.render();
            }
        },

        onNextPageClick: function () {
            this.collection.getNextPage(null);
        },

        render: function () {
            App.log('ThreadListView::render');
            this.$el.html(this.template({}));

            this.$thumbList = this.$('.thumb-list');

            return this;
        }
    });
});
