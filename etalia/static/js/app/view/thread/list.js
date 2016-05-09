define([
    'app',
    'text!app/templates/thread/list.hbs',
    'app/view/list',
    'app/view/detail',
    'app/view/ui/modal',
    'app/view/thread/detail',
    'app/view/thread/thumb',
    'app/collection/thread/thread'
], function (App, template) {

    App.View.Thread = App.View.Thread || {};

    App.View.Thread.List = App.Backbone.View.extend({
        tagName: 'div',

        template: App.Handlebars.compile(template),

        listView: null,
        detailView: null,

        events: {
            "click #thread-create-modal": "onCreateModalClick",
            "click #thread-next-page": "onNextPageClick"
        },

        initialize: function () {
            this.listenTo(this, "model:remove", this.onModelRemove);
            this.listenTo(this, "model:detail", this.openDetail);
        },

        render: function () {
            App.log('ThreadListView::render');

            this.$el.html(this.template({}));


            // Thumbs list
            this.listView = new App.View.List.create({}, {
                $target: this.$('div[data-list-placeholder]')
            });
            this.pushSubView(this.listView);

            return this;
        },

        load: function(data) {
            App.Layout.setBusy();

            this.clearList();

            if (this.collection) {
                this.collection.fullCollection.off("add", this.onCollectionAdd);
                this.collection.fullCollection.off("remove", this.onCollectionRemove);
                this.collection.fullCollection.reset();
                this.collection.reset();
            }

            this.collection = new App.Collection.Threads({
                query: data
            });

            this.collection.fullCollection.on("add", this.onCollectionAdd, this);
            this.collection.fullCollection.on("remove", this.onCollectionRemove, this);

            var that = this;
            this.collection.fetch()
                .then(function() {
                    if (that.collection.hasNextPage()) {
                        that.$('#thread-next-page').show();
                    } else {
                        that.$('#thread-next-page').hide();
                    }
                    App.Layout.setAvailable();
                });
        },

        onCollectionAdd: function (model) {
            App.log('ThreadListView::onCollectionAdd');

            // Render the thumb
            var thumbView = new App.View.Thread.Thumb({
                id: 'thread-thumb-' + model.get('id'),
                model: model,
                list: this
            });

            this.listView.addThumbView(thumbView);
        },

        onCollectionRemove: function (model) {
            App.log('ThreadListView::onCollectionRemove');

            this.listView.removeThumbById('thread-thumb-' + model.get('id'));
        },

        clearList: function() {
            this.listView.clear();
        },

        openDetail: function (model) {
            App.log('ThreadListView::onModelDetail');

            var modelIndex = this.collection.fullCollection.indexOf(model);
            if (0 > modelIndex) {
                throw 'Unexpected index';
            }

            var detailModel = new App.Model.Detail({
                list: this,
                prev: 0 < modelIndex ? this.collection.fullCollection.at(modelIndex - 1) : null,
                next: this.collection.fullCollection.at(modelIndex + 1),
                view: new App.View.Thread.Detail({
                    model: model
                })
            });

            if (this.detailView) {
                this.detailView.model.get('view').remove();
                this.detailView.model.destroy();
                this.detailView.model = detailModel;
            } else {
                this.detailView = new App.View.Detail({
                    model: detailModel
                });
                this.listenTo(this.detailView, "detail:prev", this.openDetail);
                this.listenTo(this.detailView, "detail:next", this.openDetail);
            }

            var that = this;
            model
                .fetch({data: {view: 'nested'}})
                .done(function() {
                    that.detailView.render();
                });
        },

        onModelRemove: function (model) {
            App.log('ThreadListView::onModelRemove');

            this.collection.remove(model);
            this.collection.fullCollection.remove(model);
        },

        onNextPageClick: function (e) {
            e.preventDefault();

            var that = this;
            if (this.collection.hasNextPage()) {
                this.collection.getNextPage(null).then(function() {
                    if (that.collection.hasNextPage()) {
                        that.$('#thread-next-page').show();
                    } else {
                        that.$('#thread-next-page').hide();
                    }
                });
            }
        },

        onCreateModalClick: function(e) {
            e.preventDefault();

            var that = this,
                form = App.View.Thread.CreateForm.create(),
                modal = new App.View.Ui.Modal({
                    title: 'Start a new thread',
                    content: form,
                    footer: false
                });

            form.once('validation_success', function () {
                form.model.save(null, {
                    wait: true,
                    success: function () {
                        that.collection.add(form.model, {at: 0});
                        modal.close();
                        that.openDetail(form.model);
                    },
                    error: function () {
                        // TODO
                    }
                });
            });

            form.once('cancel', function () {
                modal.close();
            });

            modal.once('hidden', function () {
                form = null;
                modal = null;
            });

            modal.render();
        }
    });


    App.View.Thread.List.create = function(options, createOptions) {
        options = options || {};

        var view = new App.View.Thread.List(options);
        if (createOptions) {
            App.View.create(view, createOptions);
        }

        return view;
    };

    return App.View.Thread.List;
});
